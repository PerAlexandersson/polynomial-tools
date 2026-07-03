//! Local weights for the Borodin--Wheeler `L`-matrix vertex model.
//!
//! This implements the face weights in Equation (3.7) of
//! Borodin--Wheeler, "Nonsymmetric Macdonald polynomials via integrable
//! vertex models".

/// A factored local face weight of the form `x^a t^b prod_m (1 - t^m)`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BorodinWheelerFaceWeight {
    pub x_power: usize,
    pub t_power: usize,
    pub one_minus_t_power_factors: Vec<usize>,
}

impl BorodinWheelerFaceWeight {
    pub fn one() -> Self {
        Self {
            x_power: 0,
            t_power: 0,
            one_minus_t_power_factors: Vec::new(),
        }
    }

    pub fn factored_string(&self) -> String {
        let mut factors = Vec::new();
        match self.x_power {
            0 => {}
            1 => factors.push("x".to_string()),
            exponent => factors.push(format!("x^{exponent}")),
        }
        for &exponent in &self.one_minus_t_power_factors {
            match exponent {
                0 => factors.push("0".to_string()),
                1 => factors.push("(1 - t)".to_string()),
                _ => factors.push(format!("(1 - t^{exponent})")),
            }
        }
        match self.t_power {
            0 => {}
            1 => factors.push("t".to_string()),
            exponent => factors.push(format!("t^{exponent}")),
        }
        if factors.is_empty() {
            "1".to_string()
        } else {
            factors.join(" ")
        }
    }
}

/// Compute the local Borodin--Wheeler `L_x(I,left;K,right)` face weight.
///
/// The color labels `left` and `right` lie in `0..=n`, while `bottom` and
/// `top` are length-`n` occupation vectors.  The result is returned in
/// factored form; `None` denotes weight zero.
pub fn borodin_wheeler_l_weight(
    bottom: &[usize],
    left: usize,
    top: &[usize],
    right: usize,
) -> Option<BorodinWheelerFaceWeight> {
    assert_eq!(
        bottom.len(),
        top.len(),
        "bottom and top occupation vectors must have the same length"
    );
    let n = bottom.len();
    assert!(
        left <= n && right <= n,
        "vertical color labels must lie in 0..=n"
    );

    if left == 0 && right == 0 {
        return (top == bottom).then(BorodinWheelerFaceWeight::one);
    }

    if left == right {
        return (top == bottom).then(|| BorodinWheelerFaceWeight {
            x_power: 1,
            t_power: suffix_sum_after(bottom, left),
            one_minus_t_power_factors: Vec::new(),
        });
    }

    if left == 0 {
        let color = right;
        return decremented_at(bottom, color).and_then(|expected_top| {
            (top == expected_top).then(|| BorodinWheelerFaceWeight {
                x_power: 1,
                t_power: suffix_sum_after(bottom, color),
                one_minus_t_power_factors: vec![bottom[color - 1]],
            })
        });
    }

    if right == 0 {
        let color = left;
        return (top == incremented_at(bottom, color)).then(BorodinWheelerFaceWeight::one);
    }

    if left < right {
        return moved_color(bottom, right, left).and_then(|expected_top| {
            (top == expected_top).then(|| BorodinWheelerFaceWeight {
                x_power: 1,
                t_power: suffix_sum_after(bottom, right),
                one_minus_t_power_factors: vec![bottom[right - 1]],
            })
        });
    }

    None
}

fn suffix_sum_after(vector: &[usize], color: usize) -> usize {
    vector[color..].iter().sum()
}

fn incremented_at(vector: &[usize], color: usize) -> Vec<usize> {
    let mut result = vector.to_vec();
    result[color - 1] += 1;
    result
}

fn decremented_at(vector: &[usize], color: usize) -> Option<Vec<usize>> {
    if vector[color - 1] == 0 {
        return None;
    }
    let mut result = vector.to_vec();
    result[color - 1] -= 1;
    Some(result)
}

fn moved_color(vector: &[usize], from_color: usize, to_color: usize) -> Option<Vec<usize>> {
    if vector[from_color - 1] == 0 {
        return None;
    }
    let mut result = vector.to_vec();
    result[from_color - 1] -= 1;
    result[to_color - 1] += 1;
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_borodin_wheeler_example_3_2() {
        let weight = borodin_wheeler_l_weight(&[1, 1, 1], 1, &[2, 0, 1], 2).unwrap();
        assert_eq!(
            weight,
            BorodinWheelerFaceWeight {
                x_power: 1,
                t_power: 1,
                one_minus_t_power_factors: vec![1],
            }
        );
        assert_eq!(weight.factored_string(), "x (1 - t) t");
    }

    #[test]
    fn matches_borodin_wheeler_example_3_3() {
        let weight = borodin_wheeler_l_weight(&[2, 1, 2], 0, &[1, 1, 2], 1).unwrap();
        assert_eq!(
            weight,
            BorodinWheelerFaceWeight {
                x_power: 1,
                t_power: 3,
                one_minus_t_power_factors: vec![2],
            }
        );
        assert_eq!(weight.factored_string(), "x (1 - t^2) t^3");
    }

    #[test]
    fn reverse_color_crossing_has_zero_weight() {
        assert_eq!(borodin_wheeler_l_weight(&[1, 1, 1], 2, &[0, 2, 1], 1), None);
    }
}
