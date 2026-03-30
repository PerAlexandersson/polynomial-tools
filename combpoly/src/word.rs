//! Words on boards: enumeration and parsing for fillings of (skew) shapes.

/// Generate all words with given content fitting on a (skew) board.
///
/// `content[j-1]` = number of copies of letter `j` (letters are 1-indexed).
/// `board[i]` = upper bound for position `i` (i.e., `w_i <= board[i]`).
/// `skew[i]` = lower bound for position `i` (i.e., `w_i > skew[i]`), defaults to 0.
pub fn words_on_board(content: &[usize], board: &[u8], skew: Option<&[u8]>) -> Vec<Vec<u8>> {
    let n = board.len();
    let k = content.len();
    let total: usize = content.iter().sum();
    assert_eq!(
        total, n,
        "Content must sum to board length (content sums to {}, board has length {})",
        total, n
    );

    let mut result = Vec::new();
    let mut current = Vec::with_capacity(n);
    let mut remaining = content.to_vec();
    let default_skew = vec![0u8; n];
    let skew = skew.unwrap_or(&default_skew);

    generate_words(
        0,
        n,
        k,
        board,
        skew,
        &mut remaining,
        &mut current,
        &mut result,
    );
    result
}

#[allow(clippy::too_many_arguments)]
fn generate_words(
    pos: usize,
    n: usize,
    k: usize,
    board: &[u8],
    skew: &[u8],
    remaining: &mut Vec<usize>,
    current: &mut Vec<u8>,
    result: &mut Vec<Vec<u8>>,
) {
    if pos == n {
        result.push(current.clone());
        return;
    }

    let lo = skew[pos] as usize + 1;
    let hi = (board[pos] as usize).min(k);

    for letter in lo..=hi {
        let j = letter - 1;
        if remaining[j] > 0 {
            remaining[j] -= 1;
            current.push(letter as u8);
            generate_words(pos + 1, n, k, board, skew, remaining, current, result);
            current.pop();
            remaining[j] += 1;
        }
    }
}

/// Parse a board string like "22233" or "2,2,2,3,3".
pub fn parse_board(s: &str) -> Vec<u8> {
    if s.contains(',') {
        s.split(',')
            .map(|x| x.trim().parse().expect("invalid number in board"))
            .collect()
    } else {
        s.chars()
            .map(|c| c.to_digit(10).expect("invalid digit in board") as u8)
            .collect()
    }
}

/// Parse a content string like "2,2,1".
pub fn parse_content(s: &str) -> Vec<usize> {
    s.split(',')
        .map(|x| x.trim().parse().expect("invalid number in content"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_words_on_board_paper_example() {
        // Content (2,2,1) on board 22233 should give 12 words
        let words = words_on_board(&[2, 2, 1], &[2, 2, 2, 3, 3], None);
        assert_eq!(words.len(), 12);
    }

    #[test]
    fn test_words_full_board() {
        // All permutations of multiset {1,1,2} on full board [2,2,2]
        let words = words_on_board(&[2, 1], &[2, 2, 2], None);
        assert_eq!(words.len(), 3); // 112, 121, 211
    }

    #[test]
    fn test_words_multinomial() {
        // All words with content (1,1,1) on board [3,3,3] = all permutations of S_3
        let words = words_on_board(&[1, 1, 1], &[3, 3, 3], None);
        assert_eq!(words.len(), 6);
    }

    #[test]
    fn test_skew_board() {
        // Content (1,1) on board [2,2], skew [1,0]: w_1 > 1 and w_2 > 0
        // So w_1 in {2} and w_2 in {1,2}, but content is one 1 and one 2
        // Only valid word: [2, 1]
        let words = words_on_board(&[1, 1], &[2, 2], Some(&[1, 0]));
        assert_eq!(words.len(), 1);
        assert_eq!(words[0], vec![2, 1]);
    }
}
