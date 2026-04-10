//! Basis enum for nonsymmetric polynomial bases.

use std::fmt;

/// The five nonsymmetric polynomial bases indexed by weak compositions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum MultiPolyBasis {
    /// Monomial basis x^α.
    Monomial,
    /// Key polynomials κ_α (Demazure characters).
    Key,
    /// Demazure atoms A_α.
    Atom,
    /// Monomial slide polynomials (Assaf-Searles 2016). Stable limit: QSym M_α.
    MonSlide,
    /// Fundamental slide polynomials (Assaf-Searles 2016). Stable limit: QSym F_α.
    FundSlide,
}

impl MultiPolyBasis {
    /// Short symbol for display.
    pub fn symbol(&self) -> &'static str {
        match self {
            MultiPolyBasis::Monomial => "x",
            MultiPolyBasis::Key => "κ",
            MultiPolyBasis::Atom => "A",
            MultiPolyBasis::MonSlide => "M",
            MultiPolyBasis::FundSlide => "F",
        }
    }
}

impl fmt::Display for MultiPolyBasis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MultiPolyBasis::Monomial => write!(f, "Monomial"),
            MultiPolyBasis::Key => write!(f, "Key"),
            MultiPolyBasis::Atom => write!(f, "Atom"),
            MultiPolyBasis::MonSlide => write!(f, "MonSlide"),
            MultiPolyBasis::FundSlide => write!(f, "FundSlide"),
        }
    }
}
